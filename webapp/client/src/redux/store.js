import { configureStore } from '@reduxjs/toolkit'
// import logger from 'redux-logger'
import coreuiReducer from './reducers/coreuiSlice'
import userReducer from './reducers/edge/userSlice'
import pageReducer from './reducers/pageSlice'
import messageReducer from './reducers/messageSlice'

// The thunk middleware was automatically added
const store = configureStore({
  reducer: {
    coreui: coreuiReducer,
    user: userReducer,
    page: pageReducer,
    message: messageReducer,
  },
  //middleware: (getDefaultMiddleware) => getDefaultMiddleware().concat(logger),
})

export default store
