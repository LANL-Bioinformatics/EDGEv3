/* eslint-disable no-unreachable */
const fs = require('fs');
const ejs = require('ejs');
const Papa = require('papaparse');
const Job = require('../edge-api/models/job');
const { nextflowConfigs, workflowList, generateWorkflowResult } = require('./workflow');
const { write2log, execCmd, sleep } = require('./common');
const logger = require('./logger');
const config = require('../config');

const generateInputs = async (projHome, projectConf, proj) => {
  // projectConf: project conf.js
  // workflowList in utils/workflow
  const workflowSettings = workflowList[projectConf.workflow.name];
  const template = String(fs.readFileSync(`${config.NEXTFLOW.TEMPLATE_DIR}/${workflowSettings.config_tmpl}`));
  const params = {
    ...projectConf.workflow.input,
    outdir: `${projHome}/${workflowSettings.outdir}`,
    projOutdir: `${projHome}/${workflowSettings.outdir}`,
    project: proj.name,
    report_config: `${config.NEXTFLOW.CONFIG_DIR}/${nextflowConfigs.report_config}`,
    slurm_config: `${config.NEXTFLOW.CONFIG_DIR}/${nextflowConfigs.slurm_config}`
  };
  if (projectConf.workflow.name === 'sra2fastq') {
    params.outdir = config.IO.SRA_BASE_DIR;
  }
  // render input template and write to nextflow_params.json
  const inputs = ejs.render(template, params);
  await fs.promises.writeFile(`${projHome}/nextflow.config`, inputs);
  return true;
};

// submit workflow - launch nextflow run
const submitWorkflow = async (proj, projectConf, inputsize) => {
  const projHome = `${config.IO.PROJECT_BASE_DIR}/${proj.code}`;
  const log = `${projHome}/log.txt`;
  // Run nextflow in <project home>/nextflow
  const workDir = `${projHome}/nextflow`;
  fs.mkdirSync(workDir);
  if (!fs.existsSync(workDir)) {
    logger.error(`Error creating directory ${workDir}:`);
    proj.status = 'failed';
    proj.updated = Date.now();
    proj.save();
    return;
  }
  // submit workflow
  const runName = `edge-${proj.code}`;
  const cmd = `cd ${workDir}; ${config.NEXTFLOW.PATH} -c ${projHome}/nextflow.config -bg -q run ${config.NEXTFLOW.WORKFLOW_DIR}/${workflowList[projectConf.workflow.name].nextflow_main} -name ${runName}`;
  write2log(log, 'Run pipeline');
  // Don't need to wait for the command to complete. It may take long time to finish and cause an error.
  // The updateJobStatus will catch the error if this command failed.
  execCmd(cmd);
  await sleep(2000); // Wait for 2 seconds
  const newJob = new Job({
    id: runName,
    project: proj.code,
    type: proj.type,
    inputsize,
    queue: 'nextflow',
    status: 'Running'
  });
  newJob.save().catch(err => { logger.error('falied to save to nextflowjob: ', err); });
  proj.status = 'running';
  proj.updated = Date.now();
  proj.save();
};

const getPid = async (proj) => {
  // To stop the running pipeline depends on the executor.
  // If is local, find pid in .nextflow.pid and kill process and all descendant processes: pkill -TERM -P <pid>
  // If is slurm, delete slurm job?
  const pidFile = `${config.IO.PROJECT_BASE_DIR}/${proj.code}/nextflow/.nextflow.pid`;
  if (fs.existsSync(pidFile)) {
    let all = fs.readFileSync(pidFile, 'utf8');
    all = all.trim();  // final crlf in file
    const lines = all.split('\n');
    if (lines[0]) {
      return parseInt(lines[0], 10);
    }
  }
  return null;
};
// check pid
const pidIsRunning = (pid) => {
  try {
    // a signal of 0 can be used to test for the existence of a process.
    process.kill(pid, 0);
    return true;
  } catch (e) {
    return false;
  }
};

const abortJob = async (proj) => {
  // To stop the running pipeline depends on the executor.
  // If is local, find pid in .nextflow.pid and kill process and all descendant processes: pkill -TERM -P <pid>
  // If is slurm, delete slurm job?
  const pid = await getPid(proj);
  if (pid && pidIsRunning(pid)) {
    const cmd = `pkill -TERM -P ${pid}`;
    // Don't need to wait for the deletion, the process may already complete
    execCmd(cmd);
  }
  // delete job
  Job.deleteOne({ project: proj.code }, (err) => {
    if (err) {
      logger.error(`Failed to delete job from DB ${proj.code}:${err}`);
    }
  });
};

const getJobMetadata = async (proj) => {
  const traceFile = `${config.IO.PROJECT_BASE_DIR}/${proj.code}/nextflow/trace.txt`;
  if (!fs.existsSync(traceFile)) {
    return [];
  }
  // get job metadata in trace.txt, convert tab delimiter file to json
  const jobMetadata = Papa.parse(fs.readFileSync(traceFile).toString(), { delimiter: '\t', header: true, skipEmptyLines: true }).data;
  return jobMetadata;
};

const generateRunStats = async (project) => {
  const stats = await getJobMetadata(project);
  fs.writeFileSync(`${config.IO.PROJECT_BASE_DIR}/${project.code}/run_stats.json`, JSON.stringify({ 'stats': stats }));
};

const getJobStatus = (statusStr) => {
  // parse output from 'nextflow log <run name> -f status
  const statuses = statusStr.split(/\n/);
  let completeCnt = 0;
  let i = 0;
  for (i = 0; i < statuses.length; i += 1) {
    const status = statuses[i].trim();
    if (status === '' || status === 'COMPLETED') {
      // empty line === COMPLETED
      completeCnt += 1;
    }
    if (status === 'ABORTED') {
      return 'Aborted';
    }
  }
  if (completeCnt === statuses.length) {
    return 'Succeeded';
  }
  return 'Failed';
};

const updateJobStatus = async (job, proj) => {
  // get job status
  const projHome = `${config.IO.PROJECT_BASE_DIR}/${proj.code}`;
  // Pipeline status. Possible values are: OK, ERR and empty
  let cmd = `cd ${projHome}/nextflow; ${config.NEXTFLOW.PATH} log|awk '/${job.id}/ &&(/OK/||/ERR/)'|awk '{split($0,array,/\t/); print array[4]}'`;
  let ret = await execCmd(cmd);

  if (!ret || ret.code !== 0) {
    // command failed
    return;
  }
  // if empty, check pid
  if (ret.message === '') {
    const pid = await getPid(proj);
    if (pid && pidIsRunning(pid)) {
      // workflow is still running, update job updated datetime to move job to the end of job queue
      job.updated = Date.now();
      job.save();
    } else {
      // workflow failed
      job.status = 'Failed';
      job.updated = Date.now();
      job.save();
      // result not as expected
      proj.status = 'failed';
      proj.updated = Date.now();
      proj.save();
      write2log(`${projHome}/log.txt`, 'Nextflow job status: Failed');
    }
    return;
  }

  // Task status. Possible values are: COMPLETED, FAILED, and ABORTED.
  cmd = `cd ${projHome}/nextflow; ${config.NEXTFLOW.PATH} log ${job.id} -f status`;
  ret = await execCmd(cmd);
  if (!ret || ret.code !== 0) {
    // command failed
    return;
  }
  // find job status
  const newStatus = getJobStatus(ret.message);
  // update project status
  if (job.status !== newStatus) {
    let status = null;
    if (newStatus === 'Aborted') {
      status = 'failed';
    } else if (newStatus === 'Succeeded') {
      // generate result.json
      logger.info('generate workflow result.json');
      try {
        generateWorkflowResult(proj);
      } catch (e) {
        job.status = newStatus;
        job.updated = Date.now();
        job.save();
        // result not as expected
        proj.status = 'failed';
        proj.updated = Date.now();
        proj.save();
        throw e;
      }
      status = 'complete';
    } else if (newStatus === 'Failed') {
      status = 'failed';
    }
    proj.status = status;
    proj.updated = Date.now();
    proj.save();
    write2log(`${projHome}/log.txt`, `Nextflow job status: ${newStatus}`);
  }
  // update job even its status unchanged. We need set new updated time for this job.
  if (newStatus === 'Aborted') {
    // delete job
    Job.deleteOne({ project: proj.code }, (err) => {
      if (err) {
        logger.error(`Failed to delete job from DB ${proj.code}:${err}`);
      }
    });
  } else {
    job.status = newStatus;
    job.updated = Date.now();
    job.save();
  }
};

module.exports = {
  generateInputs,
  submitWorkflow,
  generateRunStats,
  abortJob,
  getJobMetadata,
  updateJobStatus,
};
