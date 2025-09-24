const Upload = require('../models/upload');

const getRealName = async (dir) => {
  let name = dir;
  // check if input is a uploaded file
  const upload = await Upload.findOne({ 'code': dir });

  if (upload) {
    name = upload.name;
  }
  return name;
};

const updateUpload = async (query, req) => {
  try {
    const upload = await Upload.findOne(query);
    if (!upload) {
      return null;
    }
    upload.name = req.body.name;
    upload.desc = req.body.desc;
    upload.folder = req.body.folder;
    upload.public = req.body.public;
    upload.sharedTo = req.body.sharedTo;
    if (req.body.status) {
      upload.status = req.body.status;
    }

    upload.updated = Date.now();
    return await upload.save();
  } catch (err) {
    return Promise.reject(err);
  }
};

const getUploadedSize = async (owner) => {
  try {
    // find size of all uploaded files owned by user
    const result = await Upload.aggregate([{
      $match: { $and: [{ 'status': { $ne: 'delete' }, 'owner': { $eq: owner } }] },
    }, {
      $group: {
        _id: null,
        total: {
          $sum: '$size'
        }
      }
    }]);
    let size = 0;
    if (result[0]) {
      size = result[0].total;
    }
    return size;
  } catch (err) {
    return Promise.reject(err);
  }
};

const getUploadFolderOptions = async (owner) => {
  try {
    // get folder options
    const folders = await Upload.distinct('folder', { 'status': { $ne: 'delete' }, 'owner': { $eq: owner } });
    const options = [];
    if (folders.length > 0) {
      folders.forEach(folder => {
        options.push({ value: folder, label: folder });
      });
    } else {
      options.push({ value: 'main', label: 'main' });
    }
    return options;
  } catch (err) {
    return Promise.reject(err);
  }
};

module.exports = {
  getRealName,
  updateUpload,
  getUploadedSize,
  getUploadFolderOptions,
};
