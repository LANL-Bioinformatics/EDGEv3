const mongoose = require('mongoose');
const { uploadStatus } = require('../utils/conf');

// Create Schema
const uploadSchema = new mongoose.Schema({
  name: {
    type: String,
    required: true,
    trim: true,
  },
  desc: {
    type: String,
    trim: true,
  },
  folder: {
    type: String,
    required: true
  },
  size: {
    type: Number,
    required: true,
  },
  type: {
    type: String,
    required: true,
  },
  code: {
    type: String,
    required: true,
    unique: true,
  },
  status: {
    type: String,
    default: uploadStatus[0],
    enum: uploadStatus,
  },
  public: {
    type: Boolean,
    default: false,
  },
  owner: {
    type: String,
    required: true,
  },
  sharedTo: [{
    type: String,
  }],
}, {
  timestamps: {
    createdAt: 'created', // Use `created` to store the created date
    updatedAt: 'updated', // and `updated` to store the last updated date
  },
});

module.exports = mongoose.model('Upload', uploadSchema);
