const { exec } = require('child_process');
const moment = require('moment');

const logger = require('../utils/logger');
const config = require('../config');

module.exports = function dbBackup() {
  logger.debug('DB backup');
  // mongodump
  const dateStringWithTime = moment(new Date()).format('YYYY-MM-DD:HH:mm');
  const cmd = `mongodump --db ${config.DATABASE.NAME} --out ${config.DATABASE.BACKUP_DIR}/db-backup_${dateStringWithTime}`;

  logger.info(cmd);
  // run local
  exec(cmd, (error, stdout, stderr) => {
    if (error) {
      logger.error(error.message);
    }
    if (stderr) {
      logger.error(stderr);
    }
  });
};
