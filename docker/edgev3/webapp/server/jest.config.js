const { defaults } = require('jest-config');

module.exports = {
  // environmental variables for Jest's tests
  setupFiles: ['./tests/test.env'],
  moduleFileExtensions: [...defaults.moduleFileExtensions, 'js', 'jsx'],
};
