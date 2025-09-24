const { body, validationResult } = require('express-validator');

const validationRules = () => [
  body('name').trim().notEmpty().escape().withMessage('Invalid name.'),
  body('folder').trim().notEmpty().withMessage('Invalid folder.'),
  body('desc').optional(),
  body('type').trim().notEmpty().escape().withMessage('Invalid type'),
  body('size').trim().isInt().withMessage('Invalid size.'),
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
