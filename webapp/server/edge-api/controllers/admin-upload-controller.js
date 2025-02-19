const Upload = require('../models/upload');
const { updateUpload } = require('../utils/upload');
const logger = require('../../utils/logger');
const config = require('../../config');

const sysError = config.APP.API_ERROR;

// Update upload
const updateOne = async (req, res) => {
  try {
    logger.debug(`/api/admin/uploads/${req.params.code} update`);

    const query = { 'status': { $ne: 'delete' }, code: { $eq: req.params.code } };
    const upload = await updateUpload(query, req);

    if (!upload) {
      logger.error(`upload ${req.params.code} not found or access denied.`);
      return res.status(400).json({
        error: { upload: `upload ${req.params.code} not found or access denied` },
        message: 'Action failed',
        success: false,
      });
    }
    return res.send({
      upload,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`Admin update upload failed: ${err}`);

    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

// Find all uploads
const getAll = async (req, res) => {
  try {
    logger.debug('/api/admin/uploads');
    // find all uploads
    const uploads = await Upload.find({ 'status': { $ne: 'delete' } }).sort([['updated', -1]]);

    return res.json({
      uploads,
      message: 'Action successful',
      success: true,
    });

  } catch (err) {
    logger.error(`Admin get uploads failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

module.exports = {
  updateOne,
  getAll,
};
