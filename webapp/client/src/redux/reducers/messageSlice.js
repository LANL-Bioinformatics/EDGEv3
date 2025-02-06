import { createSlice } from '@reduxjs/toolkit'

const initialState = { errors: {}, messages: {} }

export const messageSlice = createSlice({
  name: 'message',
  initialState,
  // The `reducers` field lets us define reducers and generate associated actions
  reducers: {
    // Use the PayloadAction type to declare the contents of `action.payload`
    addMessage: (state, action) => {
      state.messages = { ...state.messages, ...action.payload }
    },
    cleanMessage: (state) => {
      state.messages = {}
    },
    addError: (state, action) => {
      // payload is in json format
      state.errors = { ...state.errors, ...action.payload }
    },
    cleanError: (state) => {
      state.errors = {}
    },
  },
})

export const { addMessage, cleanMessage, addError, cleanError } = messageSlice.actions

export default messageSlice.reducer
