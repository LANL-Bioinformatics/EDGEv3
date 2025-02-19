import { createSlice } from '@reduxjs/toolkit'

const initialState = {
  sidebarShow: true,
  sidebarUnfoldable: false,
  theme: 'light',
}

export const coreuiSlice = createSlice({
  name: 'coreui',
  initialState,
  // The `reducers` field lets us define reducers and generate associated actions
  reducers: {
    // Redux Toolkit allows us to write "mutating" logic in reducers. It
    // doesn't actually mutate the state because it uses the Immer library,
    // which detects changes to a "draft state" and produces a brand new
    // immutable state based off those changes

    setSidebar: (state, action) => {
      state.sidebarShow = action.payload
    },
    setSidebarUnfoldable: (state, action) => {
      state.sidebarUnfoldable = action.payload
    },
    setTheme: (state, action) => {
      state.theme = action.payload
    },
  },
})

export const { setSidebar, setSidebarUnfoldable, setTheme } = coreuiSlice.actions

export default coreuiSlice.reducer
