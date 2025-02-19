const bcrypt = require('bcryptjs');
const jwt = require('jsonwebtoken');
const config = require('../../config');

const encodePassword = async (password) => new Promise((resolve, reject) => {
  // Hash password before saving in database
  bcrypt.genSalt(10, (error, salt) => {
    bcrypt.hash(password, salt, (err, hash) => {
      if (err) {
        return reject(err);
      }
      return resolve(hash);
    });
  });
});

const signToken = async (payload) => new Promise((resolve, reject) => {
  // Sign token
  jwt.sign(
    payload,
    config.AUTH.JWT_SECRET,
    {
      expiresIn: 31556926, // 1 year in seconds
    },
    (err, token) => {
      if (err) {
        return reject(err);
      }
      return resolve(token);
    },
  );
});

module.exports = {
  encodePassword,
  signToken,
};
