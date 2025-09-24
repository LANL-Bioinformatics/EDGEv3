const JwtStrategy = require('passport-jwt').Strategy;
const { ExtractJwt } = require('passport-jwt');
const User = require('../models/user');
const logger = require('../../utils/logger');
const config = require('../../config');

const opts = {};
opts.jwtFromRequest = ExtractJwt.fromAuthHeaderAsBearerToken();
opts.secretOrKey = config.AUTH.JWT_SECRET;

module.exports = passport => {
  passport.use(
    'user',
    new JwtStrategy(opts, (jwtPayload, done) => {
      User.findById(jwtPayload.id)
        .then(user => {
          if (user && user.active) {
            return done(null, user);
          }
          return done(null, false);
        })
        .catch(err => {
          logger.error(err);
        });
    })
  );

  passport.use(
    'admin',
    new JwtStrategy(opts, (jwtPayload, done) => {
      User.findById(jwtPayload.id)
        .then(user => {
          if (user && user.active && user.role === 'admin') {
            return done(null, user);
          }
          return done(null, false);
        })
        .catch(err => {
          logger.error(err);
        });
    })
  );
};
