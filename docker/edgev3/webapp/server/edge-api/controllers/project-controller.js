const Project = require('../models/project');
const { getProject, getProjectConf, getProjectOutputs, getProjectResult, getProjectRunStats } = require('../utils/project');
const logger = require('../../utils/logger');
const config = require('../../config');

const sysError = config.APP.API_ERROR;

// Find all public projects
const getAll = async (req, res) => {
  try {
    logger.debug('/api/public/projects');
    const projects = await Project.find({ status: { $ne: 'delete' }, public: true }, { sharedTo: 0 }).sort([['updated', -1]]);
    return res.json({
      projects,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`/api/public/projects failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

// Find all public projects
const getOne = async (req, res) => {
  try {
    logger.debug(`/api/public/projects/${req.params.code}`);
    const project = await getProject(req.params.code, 'public');
    if (!project) {
      logger.error(`project ${req.params.code} not found or access denied.`);
      return res.status(400).json({
        error: { project: `project ${req.params.code} not found or access denied` },
        message: 'Action failed',
        success: false,
      });
    }
    return res.json({
      project,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`/api/public/projects/${req.params.code} failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

// Get project configuration file
const getConf = async (req, res) => {
  try {
    logger.debug(`/api/public/projects/${req.params.code}/conf`);
    const conf = await getProjectConf(req.params.code, 'public', req);

    return res.json({
      conf,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`/api/public/projects/${req.params.code}/conf failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

// Get project result file
const getResult = async (req, res) => {
  try {
    logger.debug(`/api/public/projects/${req.params.code}/result`);
    const result = await getProjectResult(req.params.code, 'public', req);

    return res.json({
      result,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`/api/public/projects/${req.params.code}/result failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

// Get project run stats
const getRunStats = async (req, res) => {
  try {
    logger.debug(`/api/public/projects/${req.params.code}/runStats`);
    const runStats = await getProjectRunStats(req.params.code, 'public', req);

    return res.json({
      runStats,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`/api/public/projects/${req.params.code}/runStats failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

// Get project output files
const getOutputs = async (req, res) => {
  try {
    logger.debug(`/api/public/projects/${req.params.code}/outputs`);
    const files = await getProjectOutputs(req.params.code, 'public', req);

    return res.json({
      fileData: files,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`/api/public/projects/${req.params.code}/outputs failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};
module.exports = {
  getAll,
  getOne,
  getConf,
  getOutputs,
  getResult,
  getRunStats,
};
