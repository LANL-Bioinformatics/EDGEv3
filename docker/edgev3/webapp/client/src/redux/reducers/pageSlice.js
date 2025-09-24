import { createSlice } from '@reduxjs/toolkit'

const initialState = {
  submittingForm: false,
  path: null,
}

export const pageSlice = createSlice({
  name: 'page',
  initialState,
  // The `reducers` field lets us define reducers and generate associated actions
  reducers: {
    // Use the PayloadAction type to declare the contents of `action.payload`
    setSubmittingForm: (state, action) => {
      state.submittingForm = action.payload
    },
    setPath: (state, action) => {
      state.path = action.payload
    },
  },
})

// The function below is called a selector and allows us to select a value from
// the state. Selectors can also be defined inline where they're used instead of
// in the slice file. For example: `useSelector((state: RootState) => state.counter.value)`
export const selectSubmittingForm = (state) => state.page.submittingForm

export const { setSubmittingForm, setPath } = pageSlice.actions

export default pageSlice.reducer
