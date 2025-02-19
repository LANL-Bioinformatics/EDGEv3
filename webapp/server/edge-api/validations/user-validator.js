const { body, validationResult } = require('express-validator');

const validationRules = () => [
  body('email').trim().normalizeEmail().isEmail().withMessage('Invalid email address.'),
  body('firstName').trim().matches(/^[a-zA-Z]+$/).withMessage('Invalid firstName.'),
  body('lastName').trim().matches(/^[a-zA-Z]+$/).withMessage('Invalid lastName.'),
  body('password', 'Password must be at least 8 characters long and contain at least one uppercase, at least one lower case and at least one special character.').isLength({ min: 8 })
    .trim().matches(/^(?=.*\d)(?=.*[a-z])(?=.*[A-Z])[a-zA-Z\d@$.!%*#?&]/),
  body('confirmPassword').notEmpty().withMessage('Confirm Password should not be empty.')
    .trim().custom((value, { req }) => {
      if (value !== req.body.password) {
        throw new Error('Password confirmation does not match with password.');
      }
      return true;
    }),
  body('active').optional().trim().isBoolean().withMessage('Invalid active. Must be true or false.'),
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
