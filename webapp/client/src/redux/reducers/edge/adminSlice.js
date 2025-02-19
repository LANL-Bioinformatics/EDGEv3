import { createAsyncThunk } from '@reduxjs/toolkit'

import { putData, deleteData, apis } from 'src/edge/common/util'
import { setSubmittingForm } from '../pageSlice'
import { cleanMessage, addError, cleanError } from '../messageSlice'

const updateProjectAsync = createAsyncThunk(
  'admin/updateProject',
  async (projData, { dispatch }) => {
    try {
      await putData(`${apis.adminProjects}/${projData.code}`, projData)
    } catch (err) {
      if (typeof err === 'string') {
        dispatch(addError({ [projData.code]: err }))
      } else {
        if (err.error) {
          dispatch(addError({ [projData.code]: JSON.stringify(err.error) }))
        } else {
          dispatch(addError({ [projData.code]: 'API server error' }))
        }
      }
    }
  },
)

const updateUploadAsync = createAsyncThunk(
  'admin/updateUpload',
  async (uploadData, { dispatch }) => {
    try {
      await putData(`${apis.adminUploads}/${uploadData.code}`, uploadData)
    } catch (err) {
      if (typeof err === 'string') {
        dispatch(addError({ [uploadData.code]: err }))
      } else {
        if (err.error) {
          dispatch(addError({ [uploadData.code]: JSON.stringify(err.error) }))
        } else {
          dispatch(addError({ [uploadData.code]: 'API server error' }))
        }
      }
    }
  },
)

const updateUserAsync = createAsyncThunk('admin/updateUser', async (userData, { dispatch }) => {
  try {
    await putData(`${apis.adminUsers}/${userData.email}`, userData)
    dispatch(setSubmittingForm(false))
  } catch (err) {
    if (typeof err === 'string') {
      dispatch(addError({ update: err }))
    } else {
      if (err.error) {
        dispatch(addError({ update: JSON.stringify(err.error) }))
      } else {
        dispatch(addError({ update: 'API server error' }))
      }
    }
    dispatch(setSubmittingForm(false))
  }
})

const deleteUserAsync = createAsyncThunk('admin/deleteUser', async (userData, { dispatch }) => {
  try {
    await deleteData(`${apis.adminUsers}/${userData.email}`, userData)
    dispatch(setSubmittingForm(false))
  } catch (err) {
    if (typeof err === 'string') {
      dispatch(addError({ delete: err }))
    } else {
      if (err.error) {
        dispatch(addError({ delete: JSON.stringify(err.error) }))
      } else {
        dispatch(addError({ delete: 'API server error' }))
      }
    }
    dispatch(setSubmittingForm(false))
  }
})

export const updateProjectAdmin = (projData) => (dispatch, getState) => {
  dispatch(updateProjectAsync(projData))
}

export const updateUploadAdmin = (uploadData) => (dispatch, getState) => {
  dispatch(updateUploadAsync(uploadData))
}

export const updateUserAdmin = (userData) => (dispatch, getState) => {
  dispatch(cleanMessage())
  dispatch(cleanError())
  dispatch(setSubmittingForm(true))
  dispatch(updateUserAsync(userData))
}

export const deleteUserAdmin = (userData) => (dispatch, getState) => {
  dispatch(cleanMessage())
  dispatch(cleanError())
  dispatch(setSubmittingForm(true))
  dispatch(deleteUserAsync(userData))
}
