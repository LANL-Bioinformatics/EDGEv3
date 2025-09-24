const fs = require('fs');
const moment = require('moment');
const Project = require('../edge-api/models/project');
const Job = require('../edge-api/models/job');
const logger = require('../utils/logger');
const config = require('../config');

module.exports = async function projectDeletionMonitor() {
  logger.debug('project deletion monitor');
  try {
    // delete file after deleteGracePeriod
    const deleteGracePeriod = moment().subtract(config.CRON.PROJECT_DELETION_GRACE_PERIOD_DAYS, 'days');
    const projs = await Project.find({ 'status': 'delete', 'updated': { '$lte': deleteGracePeriod } });
    let i;
    for (i = 0; i < projs.length; i += 1) {
      const { code } = projs[i];
      logger.info(`delete project: ${code}`);
      const path = `${config.IO.PROJECT_BASE_DIR}/${code}`;
      // delete directory recursively
      try {
        fs.rmdirSync(path, { recursive: true });
        logger.info(`deleted ${path}`);
      } catch (err) {
        logger.error(`Failed to delete ${path}:${err}`);
        return;
      }
      // delete from database
      Project.deleteOne({ code }, (err) => {
        if (err) {
          logger.error(`Failed to delete project from DB ${code}:${err}`);
        }
      });
      Job.deleteOne({ project: code }, (err) => {
        if (err) {
          logger.error(`Failed to delete job from DB ${code}:${err}`);
        }
      });
    }
  } catch (err) {
    logger.error(`projectMonitor failed:${err}`);
  }
};
