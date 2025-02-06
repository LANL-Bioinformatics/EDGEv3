const Project = require('../edge-api/models/project');
const Job = require('../edge-api/models/job');
const { abortJob, updateJobStatus } = require('../utils/nextflow');
const logger = require('../utils/logger');

module.exports = async function nextflowJobMonitor() {
  logger.debug('nextflow job monitor');
  try {
    // only process one job at each time based on job updated time
    const jobs = await Job.find({ 'queue': 'nextflow', 'status': { $in: ['Submitted', 'Running'] } }).sort({ updated: 1 });
    const job = jobs[0];
    if (!job) {
      logger.debug('No nextflow job to process');
      return;
    }
    logger.debug(`nextflow ${job.id}`);
    // find related project
    const proj = await Project.findOne({ 'code': job.project });
    if (proj) {
      if (proj.status === 'delete') {
        // abort job
        abortJob(proj);
      } else {
        await updateJobStatus(job, proj);
      }
    } else {
      // delete from database
      Job.deleteOne({ project: job.project }, (err) => {
        if (err) {
          logger.error(`Failed to delete job from DB ${job.project}:${err}`);
        }
      });
    }
  } catch (err) {
    logger.error(`nextflowJobMonitor failed:${err}`);
  }
};

