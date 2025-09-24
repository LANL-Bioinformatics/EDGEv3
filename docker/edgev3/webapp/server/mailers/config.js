const nodemailer = require('nodemailer');
const mg = require('nodemailer-mailgun-transport');
const Email = require('email-templates');
const config = require('../config');

const mailConfig = {};
let transporter = null;
if (config.EMAIL.SERVICE_PROXY_URL) {
  mailConfig.proxy = config.EMAIL.SERVICE_PROXY_URL;
}
mailConfig.auth = {
  api_key: config.EMAIL.MAILGUN_API_KEY,
  domain: config.EMAIL.MAILGUN_DOMAIN,
};
transporter = nodemailer.createTransport(mg(mailConfig));

const email = new Email({
  views: {
    root: 'email_templates',
    options: { extension: 'ejs' },
  },
  message: {
    from: config.EMAIL.FROM_ADDRESS,
  },
  send: true,
  transport: transporter,
});

// Variables in templates
const vars = {
  'getActivationLink': {
    subject: config.EMAIL.ACTIVATE_USER_SUBJECT,
    actionName: config.EMAIL.ACTIVATE_USER_ACTION,
    actionMessage: config.EMAIL.ACTIVATE_USER_ACTION_MESSAGE,
  },
  'getResetPasswordLink': {
    subject: config.EMAIL.RESETPASSWORD_SUBJECT,
    actionName: config.EMAIL.RESETPASSWORD_ACTION,
    actionMessage: config.EMAIL.RESETPASSWORD_ACTION_MESSAGE,
  },
  'projectStatusNotification': {
    subject: config.EMAIL.PROJECT_STATUS_SUBJECT,
  }
};

module.exports = {
  email,
  vars,
};
