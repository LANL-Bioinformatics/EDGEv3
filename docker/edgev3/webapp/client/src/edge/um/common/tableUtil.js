import { createTheme } from '@mui/material'
import { colors } from '../../common/util'

export const theme = createTheme({
  colorSchemes: { light: true, dark: true },
  palette: {
    primary: {
      main: colors.primary,
      dark: '#FFFFFF',
      light: '#000000',
      contrastText: '#fff',
    },
    secondary: {
      main: colors.info,
    },
  },
})

export const projectStatusColors = {
  running: 'warning',
  'in queue': 'secondary',
  failed: 'danger',
  complete: 'success',
  rerun: 'info',
  processing: 'primary',
  submitted: 'primary',
}

export const projectStatusNames = {
  running: 'Running',
  'in queue': 'In queue',
  failed: 'Failed',
  complete: 'Complete',
  rerun: 'Re-run',
  processing: 'Processing',
  submitted: 'Submitted',
}

export const userStatusColors = {
  false: 'secondary',
  true: 'success',
}

export const userStatusNames = {
  false: 'Inactive',
  true: 'Active',
}

export const userTypeColors = {
  user: 'primary',
  admin: 'danger',
}

export const userTypeNames = {
  user: 'User',
  admin: 'Admin',
}

const strongPassword = new RegExp('(?=.*[a-z])(?=.*[A-Z])(?=.*[0-9])(?=.*[^A-Za-z0-9])(?=.{8,})')
export const validatePassword = (value) => !value.length || strongPassword.test(value)
export const validateRole = (value) => ['admin', 'user'].includes(value)
export const validateRequired = (value) => !!value.length
export const validateBoolean = (value) => ['true', 'false'].includes(value)
