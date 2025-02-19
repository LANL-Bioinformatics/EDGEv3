const axios = require('axios');
const ufs = require('url-file-size');
const fs = require('fs');
const path = require('path');
const { exec } = require('child_process');
const logger = require('./logger');

// append message to a log
const write2log = (log, msg) => {
  fs.appendFile(log, `${msg}\n`, (err) => {
    if (err) logger.error(`Failed to write to ${log}: ${err}`);
  });
};

// post data
const postData = (url, params, header) => new Promise((resolve, reject) => {
  axios
    .post(url, params, header)
    .then(response => {
      const { data } = response;
      resolve(data);
    })
    .catch(err => {
      if (err.response) {
        reject(err.response);
      } else {
        reject(err);
      }
    });

});

// get data
const getData = (url) => new Promise((resolve, reject) => {
  axios
    .get(url)
    .then(response => {
      const { data } = response;
      resolve(data);
    })
    .catch(err => {
      if (err.response) {
        reject(err.response);
      } else {
        reject(err);
      }
    });

});

const timeFormat = (d) => {
  let hours = d.days() * 24 + d.hours();
  if (hours < 10) {
    hours = `0${hours}`;
  }
  let minutes = d.minutes();
  if (minutes < 10) {
    minutes = `0${minutes}`;
  }
  let seconds = d.seconds();
  if (seconds < 10) {
    seconds = `0${seconds}`;
  }
  return `${hours}:${minutes}:${seconds}`;
};

// Get file data for file browser
// Return all files matched the extentions in a directory and sub directories
const getAllFiles = (dirPath, arrayOfFilesIn, extentions, displayPath, apiPath, fileRelPath) => {
  const files = fs.readdirSync(dirPath);

  let arrayOfFiles = arrayOfFilesIn || [];

  files.forEach((file) => {
    try {
      if (fs.statSync(`${dirPath}/${file}`).isDirectory()) {
        arrayOfFiles = getAllFiles(`${dirPath}/${file}`, arrayOfFiles, extentions, `${displayPath}/${file}`, `${apiPath}/${file}`, `${fileRelPath}/${file}`);
      } else {
        let pass = false;
        if (extentions && extentions.length > 0) {
          for (let i = 0; i < extentions.length; i += 1) {
            if (file === extentions[i]) {
              pass = true;
            } else if (file.endsWith(`.${extentions[i]}`)) {
              pass = true;
            }
          }
        } else {
          // get all files
          pass = true;
        }
        if (pass) {
          const newFilePath = path.join(dirPath, '/', file);
          const newFileRelPath = path.join(fileRelPath, '/', file);
          const stats = fs.statSync(newFilePath);
          const newDisplayPath = path.join(displayPath, '/', file);
          const newApiPath = path.join(apiPath, '/', file);
          arrayOfFiles.push({
            key: newDisplayPath,
            name: file,
            path: newApiPath,
            url: newApiPath,
            filePath: newFileRelPath,
            size: stats.size,
            modified: Number(new Date(stats.mtime)),
          });
        }
      }
    } catch (err) {
      logger.error(`getAllFile failed: ${err}`);
    }
  });

  return arrayOfFiles;
};

const fileStats = async (file) => {
  let stats = {};
  if (!file) {
    return { size: 0 };
  }
  if (file.toLowerCase().startsWith('http')) {
    stats = await ufs(file)
      .then(size => ({ size }))
      .catch(() => ({ size: 0 }));
  } else {
    stats = fs.statSync(file);
  }
  return stats;
};

const findInputsize = async (projectConf) => {
  if (!projectConf.files) {
    return 0;
  }
  let size = 0;
  await Promise.all(projectConf.files.map(async (file) => {
    if (file !== '') {
      // not optional file without input
      const stats = await fileStats(file);
      size += stats.size;
    }
  }));
  // console.log('file size', size);
  return size;
};

const execCmd = (cmd) => new Promise((resolve, reject) => {
  // run local
  logger.info(`exec: ${cmd}`);
  exec(cmd, (error, stdout, stderr) => {
    if (error) {
      logger.error(error);
      reject(error);
    }
    if (stderr) {
      logger.error(`exec stderr: ${stderr}`);
      resolve({ code: -1, message: stderr });
    }
    logger.info(`exec stdout: ${stdout}`);
    resolve({ code: 0, message: stdout });
  });
});

const sleep = (ms) => new Promise(resolve => { setTimeout(resolve, ms); });

module.exports = {
  write2log,
  postData,
  getData,
  timeFormat,
  getAllFiles,
  findInputsize,
  execCmd,
  sleep,
};
