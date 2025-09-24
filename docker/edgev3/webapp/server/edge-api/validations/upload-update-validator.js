const { body, param, validationResult } = require('express-validator');
const { uploadStatus } = require('../utils/conf');

const validationRules = () => [
  param('code').trim().isLength({ min: 16 }).withMessage('Invalid upload code.'),
  body('name').optional().trim().notEmpty().escape().withMessage('Invalid name.'),
  body('desc').optional().trim(),
  body('folder').optional().trim().notEmpty().withMessage('Invalid folder.'),
  body('status').optional().trim().isIn(uploadStatus).withMessage('Invalid status.'),
  body('public').optional().trim().isBoolean().withMessage('Invalid public. Must be true or false.'),
  body('sharedTo').optional().isArray().withMessage('Invalid sharedTo. Must be an array of emails.'),
  body('sharedTo.*').optional().trim().normalizeEmail().isEmail().withMessage('Invalid email address'),
];

const validate = (req, res, next) => {
  const errors = validationResult(req);
  if (errors.isEmpty()) {
    return next();
  }

  const resultErrors = { error: {}, message: 'Validation failed', success: false };
  errors.array().forEach(err => { resultErrors.error[err.param] = err.msg; });
  return res.status(400).json(resultErrors);
};

module.exports = {
  validationRules,
  validate,
};
