const fs = require('fs');
const Project = require('../edge-api/models/project');
const Job = require('../edge-api/models/job');
const common = require('../utils/common');
const logger = require('../utils/logger');
const { nextflowWorkflows, workflowList } = require('../utils/workflow');
const { generateInputs, submitWorkflow } = require('../utils/nextflow');

const config = require('../config');

module.exports = async function nextflowWorkflowMonitor() {
  logger.debug('Nextflow workflow monitor');
  try {
    // only process one job at each time based on job updated time
    const jobs = await Job.find({ 'queue': 'nextflow', 'status': { $in: ['Submitted', 'Running'] } }).sort({ updated: 1 });
    // submit request only when the current nextflow running jobs less than the max allowed jobs
    if (jobs.length >= config.NEXTFLOW.NUM_JOBS_MAX) {
      return;
    }
    // get current running/submitted projects' input size
    let jobInputsize = 0;
    jobs.forEach(job => {
      jobInputsize += job.inputSize;
    });
    // only process one request at each time
    const projs = await Project.find({ 'type': { $in: nextflowWorkflows }, 'status': 'in queue' }).sort({ updated: 1 });
    const proj = projs[0];
    if (!proj) {
      logger.debug('No nextflow workflow request to process');
      return;
    }
    // parse conf.json
    const projHome = `${config.IO.PROJECT_BASE_DIR}/${proj.code}`;
    const projectConf = JSON.parse(fs.readFileSync(`${projHome}/conf.json`));

    // check input size
    const inputsize = await common.findInputsize(projectConf);
    if (inputsize > config.NEXTFLOW.JOBS_INPUT_MAX_SIZE_BYTES) {
      logger.debug(`Project ${proj.code} input size exceeded the limit.`);
      // fail project
      proj.status = 'failed';
      proj.updated = Date.now();
      proj.save();
      common.write2log(`${config.IO.PROJECT_BASE_DIR}/${proj.code}/log.txt`, 'input size exceeded the limit.');
      return;
    }
    if ((jobInputsize + inputsize) > config.NEXTFLOW.JOBS_INPUT_MAX_SIZE_BYTES) {
      logger.debug('Nextflow is busy.');
      return;
    }

    logger.info(`Processing workflow request: ${proj.code}`);
    // set project status to 'processing'
    proj.status = 'processing';
    proj.updated = Date.now();
    await proj.save();
    common.write2log(`${config.IO.PROJECT_BASE_DIR}/${proj.code}/log.txt`, 'Generate params json');
    logger.info('Generate params json');
    // process request
    // create output directory
    fs.mkdirSync(`${projHome}/${workflowList[projectConf.workflow.name].outdir}`, { recursive: true });
    // in case nextflow needs permission to write to the output directory
    fs.chmodSync(`${projHome}/${workflowList[projectConf.workflow.name].outdir}`, '777');
    // generate params.json
    await generateInputs(projHome, projectConf, proj);
    // submit workflow to nextflow
    common.write2log(`${config.IO.PROJECT_BASE_DIR}/${proj.code}/log.txt`, 'submit workflow to nextflow');
    logger.info('submit workflow to nextflow');
    submitWorkflow(proj, projectConf, inputsize);
    logger.info('Done workflow submission');
  } catch (err) {
    logger.error(`nextflowWorkflowMonitor failed:${err}`);
  }
};
