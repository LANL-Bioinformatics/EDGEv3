const Project = require('../edge-api/models/project');
const Job = require('../edge-api/models/job');
const { abortJob, updateJobStatus } = require('../utils/cromwell');
const logger = require('../utils/logger');

module.exports = async function cromWellJobMonitor() {
  logger.debug('cromwell job monitor');
  try {
    // only process one job at each time based on job updated time
    const jobs = await Job.find({ 'queue': 'cromwell', 'status': { $in: ['Submitted', 'Running'] } }).sort({ updated: 1 });

    const job = jobs[0];
    if (!job) {
      logger.debug('No cromwell job to process');
      return;
    }
    logger.debug(`cromwell ${job}`);
    // find related project
    const proj = await Project.findOne({ 'code': job.project });
    if (proj) {
      if (proj.status === 'delete') {
        // abort job
        abortJob(job);
      } else {
        updateJobStatus(job, proj);
      }
    } else {
      // delete from database
      await Job.deleteOne({ project: job.project });
    }
  } catch (err) {
    logger.error(`cromwellJobMonitor failed:${err}`);
  }
};
