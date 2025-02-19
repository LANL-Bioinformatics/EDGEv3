const { body, validationResult } = require('express-validator');

const validationRules = () => [
  body('firstName').trim().matches(/^[a-zA-Z]+$/).withMessage('Invalid firstName.'),
  body('lastName').trim().matches(/^[a-zA-Z]+$/).withMessage('Invalid lastName.'),
  body('email').trim().normalizeEmail().isEmail().withMessage('Invalid email address.'),
  body('oauth').trim().notEmpty().withMessage('Oauth should not be empty.'),
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
