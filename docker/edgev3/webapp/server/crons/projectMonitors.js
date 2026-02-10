const fs = require('fs');
const moment = require('moment');
const Project = require('../edge-api/models/project');
const Job = require('../edge-api/models/job');
const User = require('../edge-api/models/user');
const logger = require('../utils/logger');
const { projectStatusSender } = require('../mailers/senders');
const config = require('../config');

const projectMonitor = async () => {
  logger.debug('project status monitor');
  try {
    // notify complete/failed projects
    const projs = await Project.find({ 'notified': false });
    projs.forEach(proj => {
      if (config.EMAIL.SEND_PROJECT_STATUS_EMAILS) {
        if (proj.status === 'complete' || proj.status === 'failed') {
          User.findOne({ email: proj.owner }).then(user => {
            if (!user) {
              logger.debug(`User not found: ${proj.owner}`);
            } else if (user.notification.isOn) {
              const data = {
                projectName: proj.name,
                projectDesc: proj.desc,
                projectCreated: proj.created,
                projectType: proj.label,
                projectStatus: proj.status,
                projectPageURL: `${config.APP.UI_BASE_URL}/user/project?code=${proj.code}`
              };
              // logger.debug(`Notify user: ${user.notification.email}`);
              projectStatusSender(user.notification.email, data);
            }
          }).catch((err) => {
            logger.error(err);
          });
          proj.notified = true;
          proj.updated = Date.now();
          proj.save();
        }
      } else {
        proj.notified = true;
        proj.updated = Date.now();
        proj.save();
      }
    });
  } catch (err) {
    logger.error(`projectStatusMonitor failed:${err}`);
  }
};


const projectDeletionMonitor = async () => {
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

module.exports = {
  projectMonitor,
  projectDeletionMonitor,
};
