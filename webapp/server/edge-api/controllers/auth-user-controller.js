const randomize = require('randomatic');
const User = require('../models/user');
const { encodePassword, signToken } = require('../utils/user');
const Project = require('../models/project');
const logger = require('../../utils/logger');
const config = require('../../config');

const sysError = config.APP.API_ERROR;

// Send out info for user
const info = async (req, res) => {
  try {
    logger.debug(`/api/auth-user/info: ${JSON.stringify(req.user.email)}`);
    if (config.APP.SYSTEM_MESSAGE) {
      return res.send({
        info: {
          allowNewRuns: false,
          message: config.APP.SYSTEM_MESSAGE,
        },
        message: 'Action successful',
        success: true,
      });
    }
    // Get all projects with status in ('in queue', 'running', 'submitted')
    const projects = await Project.find({ 'owner': req.user.email, 'status': { $in: ['in queue', 'running', 'processing', 'submitted'] } });
    if (projects.length >= req.user.job.limit) {
      return res.send({
        info: {
          allowNewRuns: false,
          message: 'You have reached the maximum number of running projects allowed. Please wait for running projects to complete before submitting a new project.',
        },
        message: 'Action successful',
        success: true,
      });
    }
    return res.send({
      info: {
        allowNewRuns: true,
        jobLimit: req.user.job.limit - projects.length,
      },
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`System check failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

// Update user
const updateOne = async (req, res) => {
  try {
    const data = req.body;
    logger.debug(`/api/auth-user update: ${JSON.stringify(req.user.email)}`);
    // check if the email alreay exists
    const user = await User.findOne({ email: req.user.email });
    if (!user) {
      logger.debug(`Update user failed: Email ${req.user.email} not found`);
      return res.status(400).json({
        error: { update: `Email ${req.user.email} not found` },
        message: 'Action failed',
        success: false,
      });
    }
    // encode password
    let hashedPassword = null;
    if (data.password) {
      hashedPassword = await encodePassword(data.password);
    }
    // change code
    user.code = randomize('0', 6);
    // update fields
    user.firstName = req.body.firstName;
    user.lastName = req.body.lastName;
    user.notification.isOn = req.body.notification;

    if (req.body.active) {
      user.active = req.body.active;
    }
    if (req.body.role) {
      user.role = req.body.role;
    }
    if (req.body.notificationEmail) {
      user.notification.email = req.body.notificationEmail;
    }
    if (req.body.jobLimit) {
      user.job.limit = req.body.jobLimit;
    }
    if (req.body.jobPriority) {
      user.job.priority = req.body.jobPriority;
    }

    user.updated = Date.now();
    if (hashedPassword) {
      user.password = hashedPassword;
    }

    await user.save();
    logger.info('Update user successful');

    // Create JWT Payload
    const payload = {
      id: user.id,
      firstName: user.firstName,
      lastName: user.lastName,
      email: user.email,
      role: user.role,
      active: user.active,
      code: user.code,
      notification: user.notification,
    };
    // Sign token
    const token = await signToken(payload);
    return res.json({
      token: `Bearer ${token}`,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`Update user failed: ${err}`);
    // Delete the new user if exists
    await User.deleteOne({ email: req.body.email });

    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

// Find all active users for share/unshare actions
const getUsers = async (req, res) => {
  try {
    logger.debug(`/api/auth-user/users: ${JSON.stringify(req.user.email)}`);
    const users = await User.find({ active: true, email: { $ne: req.user.email } }, { firstName: 1, lastName: 1, email: 1, _id: 0 }).sort([['firstName', 1], ['lastName', 1]]);
    return res.send({
      users,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`Get users failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

module.exports = {
  info,
  updateOne,
  getUsers,
};
