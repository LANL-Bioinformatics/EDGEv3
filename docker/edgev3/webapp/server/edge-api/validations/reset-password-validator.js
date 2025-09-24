const { body, validationResult } = require('express-validator');

const validationRules = () => [
  body('email').trim().normalizeEmail().isEmail().withMessage('Invalid email address.'),
  body('token').trim().notEmpty().withMessage('Token should not be empty.'),
  body('newPassword', 'Password must be at least 8 characters long and contain at least one uppercase, at least one lower case and at least one special character.').isLength({ min: 8 })
    .trim().matches(/^(?=.*\d)(?=.*[a-z])(?=.*[A-Z])[a-zA-Z\d@$.!%*#?&]/),
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
