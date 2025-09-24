const Project = require('../edge-api/models/project');
const User = require('../edge-api/models/user');
const logger = require('../utils/logger');
const { projectStatusSender } = require('../mailers/senders');
const config = require('../config');

module.exports = async function projectMonitor() {
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
