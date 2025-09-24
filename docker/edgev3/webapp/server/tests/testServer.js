require('dotenv').config();
const express = require('express');
const cors = require('cors');
const helmet = require('helmet');
const passport = require('passport');
const indexRouter = require('../edge-api/routes/index');

const app = express();
app.use(cors({ origin: '*' }));
// Helmet helps to secure Express apps by setting various HTTP headers.
app.use(helmet());
app.use(express.urlencoded({ extended: false }));
app.use(express.json());

// Passport middleware
app.use(passport.initialize());
// Passport config
require('../edge-api/utils/passport')(passport);
// APIs
app.use('/api', indexRouter);

module.exports = app;
