const fs = require('fs');
const logger = require('../../utils/logger');
const config = require('../../config');

const sysError = config.APP.API_ERROR;

// Find all public projects
const getReflist = async (req, res) => {
  try {
    logger.debug('/api/workflow/metag/reflist');
    const rawdata = fs.readFileSync(config.WORKFLOW.REF_LIST);
    const reflist = JSON.parse(rawdata).tree;

    return res.json({
      reflist,
      message: 'Action successful',
      success: true,
    });
  } catch (err) {
    logger.error(`/api/workflow/metag/reflist failed: ${err}`);
    return res.status(500).json({
      message: sysError,
      success: false,
    });
  }
};

module.exports = {
  getReflist,
};
