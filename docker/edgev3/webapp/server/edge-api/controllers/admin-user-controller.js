const randomize = require('randomatic');
const User = require('../models/user');
const Project = require('../models/project');
const Upload = require('../models/upload');
const { encodePassword } = require('../utils/user');
const logger = require('../../utils/logger');
const config = require('../../config');

const sysError = config.APP.API_ERROR;

// Insert new user
const addOne = async (req, res) => {
  try {
    const data = req.body;
    logger.debug(`/api/admin/users add: ${JSON.stringify(data.email)}`);
    // check if the email alreay exists
    const user = await User.findOne({ email: { $eq: data.email } });
    if (user) {
      logger.debug(`Admin add user failed: Email ${data.email} already exists`);
      return res.status(400).json({
        error: { email: `Email ${data.email} already exists` },
        message: 'Action failed',
        success: false,
      });
    }
    // encode password and add new user to DB
    const hashedPassword = await encodePassword(data.password);
    const code = randomize('0', 6);
    const newUser = new User({
      ...data,
      password: hashedPassword,
      code,
      notification: { email: data.email },
    });

    await newUser.save();
    logger.info('User added successful');
    // success
    const result = {
      message: 'Action successful',
      success: true,
    };
    // return user for other api tests
    if (data.test) {
      result.user = newUser;
    }
    return res.json(result);
  } catch (err) {
    logger.error(`Admin add user failed: ${err}`);
    // Delete the new user if exists
    await User.deleteOne({ email: req.body.email });

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
    logger.debug(`/api/admin/users/${req.params.user} update`);
    // check if the email alreay exists
    const user = await User.findOne({ email: req.params.user });
    if (!user) {
      logger.debug(`Admin update user failed: Email ${req.params.user} not found`);
      return res.status(400).json({
        error: { update: `Email ${req.params.user} not found` },
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
    if (req.body.active) {
      user.active = req.body.active;
    }
    if (req.body.role) {
      user.role = req.body.role;
    }
    if (req.body.notificationEmail) {
      user.notification.email = req.body.notificationEmail;
    }
    if (req.body.notification) {
      user.notification.isOn = req.body.notification;
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
    logger.info('Admin update user successful');

    return res.json({
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`Admin update user failed: ${err}`);
    // Delete the new user if exists
    await User.deleteOne({ email: req.body.email });

    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

// Delete user
const removeOne = async (req, res) => {
  try {
    logger.debug(`/api/admin/users/${req.params.user} delete`);
    // check if the email alreay exists
    const user = await User.findOne({ email: req.params.user });
    if (!user) {
      logger.debug(`Admin delete user failed: Email ${req.params.user} not found`);
      return res.status(400).json({
        error: { delete: `Email ${req.params.user} not found` },
        message: 'Action failed',
        success: false,
      });
    }
    // find all projects associated with this user
    const projects = await Project.find({ 'status': { $ne: 'delete' }, owner: user.email }).count();
    if (projects === 1) {
      return res.status(400).json({
        error: { delete: `There is 1 project associated with user ${req.params.user}` },
        message: 'Action failed',
        success: false,
      });
    } if (projects > 1) {
      return res.status(400).json({
        error: { delete: `There are ${projects} projects associated with user ${req.params.user}` },
        message: 'Action failed',
        success: false,
      });
    }
    // check uploaded files
    const uploads = await Upload.find({ 'status': { $ne: 'delete' }, owner: user.email }).count();
    if (uploads === 1) {
      return res.status(400).json({
        delete: `There is 1 upload associated with user ${req.params.user}`,
        message: 'Action failed',
        success: false,
      });
    } if (uploads > 1) {
      return res.status(400).json({
        error: { delete: `There are ${uploads} uploads associated with user ${req.params.user}` },
        message: 'Action failed',
        success: false,
      });
    }
    // delete user
    await User.deleteOne({ email: user.email });

    return res.json({
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`Admin delete user failed: ${err}`);
    // Delete the new user if exists
    await User.deleteOne({ email: req.body.email });

    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

// Find all users
const getAll = async (req, res) => {
  try {
    logger.debug('/auth/admin/users');
    const users = await User.find({}).sort([['updated', -1]]);
    return res.send({
      users,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`Admin get users failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

module.exports = {
  addOne,
  updateOne,
  removeOne,
  getAll,
};
