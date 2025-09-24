const fs = require('fs');
const moment = require('moment');
const Upload = require('../edge-api/models/upload');
const logger = require('../utils/logger');
const config = require('../config');

module.exports = async function fileUploadMonitor() {
  logger.debug('file upload monitor');
  try {
    // delete file after deleteGracePeriod
    const deleteGracePeriod = moment().subtract(config.FILE_UPLOADS.DELETION_GRACE_PERIOD_DAYS, 'days');
    let uploads = await Upload.find({ 'status': 'delete', 'updated': { '$lte': deleteGracePeriod } });
    let i;
    for (i = 0; i < uploads.length; i += 1) {
      const { code } = uploads[i];
      // delete file
      const path = `${config.IO.UPLOADED_FILES_DIR}/${code}`;
      fs.unlink(path, (err) => {
        if (err) {
          logger.error(`Failed to delete ${path}:${err}`);
          return;
        }

        logger.info(`deleted ${path}`);
      });
      // delete from database
      Upload.deleteOne({ code }, (err) => {
        if (err) {
          logger.info(`Failed to delete upload ${code}:${err}`);
        }
      });
    }

    // change status to 'delete' if upload is older than daysKept
    const daysKept = moment().subtract(config.FILE_UPLOADS.FILE_LIFETIME_DAYS, 'days');
    uploads = await Upload.find({ status: 'live', 'created': { '$lte': daysKept } });
    for (i = 0; i < uploads.length; i += 1) {
      const { code } = uploads[i];
      logger.info(`mark file for deleting ${code}`);
      Upload.findOne({ code }).then(upload => {
        if (!upload) {
          logger.error(`Upload not found ${code}`);
        } else {
          upload.status = 'delete';
          upload.updated = Date.now();
          upload.save().then().catch(err => {
            logger.error(`Failed to update upload ${code}:${err}`);
          });
        }
      });
    }
  } catch (err) {
    logger.error(`fileUploadMonitor failed:${err}`);
  }
};
