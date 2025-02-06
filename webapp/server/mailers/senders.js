const { email, vars } = require('./config');
const logger = require('../utils/logger');

const userActionSender = (recipient, action, userName, link) => {
  const { subject } = vars[action];
  const { actionName } = vars[action];
  const { actionMessage } = vars[action];

  email
    .send({
      template: 'user_action',
      message: {
        to: recipient,
      },
      locals: {
        subject,
        userName,
        actionName,
        link,
        actionMessage,
      },
    })
    .then((res) => {
      logger.info(res);
    })
    .catch((err) => { logger.error(`${action} mail sent to ${recipient}: failed: ${err}`); });
};

const projectStatusSender = (recipient, data) => {
  const action = 'projectStatusNotification';
  const { subject } = vars[action];

  email
    .send({
      template: 'project_notification',
      message: {
        to: recipient,
      },
      locals: {
        subject, ...data
      },
    })
    .then((res) => {
      logger.info(res);
    })
    .catch((err) => { logger.error(`${action} mail sent to ${recipient}: failed: ${err}`); });
};

module.exports = {
  userActionSender,
  projectStatusSender
};
