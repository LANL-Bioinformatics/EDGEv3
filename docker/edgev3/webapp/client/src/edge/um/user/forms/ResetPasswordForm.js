import React from 'react'
import { Form, Input, Button, InputGroup, InputGroupText } from 'reactstrap'
import { useForm } from 'react-hook-form'
import CIcon from '@coreui/icons-react'
import { cilLockLocked } from '@coreui/icons'
import { MyTooltip } from '../../../common/MyTooltip'
import { umTips } from '../../defaults'

const ResetPasswordForm = (props) => {
  const {
    register,
    handleSubmit,
    formState: { errors },
    watch,
  } = useForm({
    mode: 'onChange',
  })

  const passwordReg = {
    ...register('newpassword', {
      required: 'Password is required',
      minLength: { value: 8, message: 'Must be at least 8 characters long' },
      validate: {
        hasUpperCase: (value) =>
          /[A-Z]/.test(value) || 'Must contain at lease one uppercase letter',
        hasLowerCase: (value) =>
          /[a-z]/.test(value) || 'Must contain at lease one lowercase letter',
        hasNumber: (value) => /[0-9]/.test(value) || 'Must contain at lease one number',
        hasSpecialChar: (value) =>
          /[^A-Za-z0-9 ]/.test(value) || 'Must contain at lease one special character',
      },
    }),
  }

  const password2Reg = {
    ...register('newpassword2', {
      validate: (value) => value === watch('newpassword') || 'The passwords do not match',
    }),
  }

  return (
    <Form onSubmit={handleSubmit(props.onSubmit)}>
      <span className="red-text">{props.errors.resetPassword}</span>
      <span className="red-text">{props.messages.resetPassword}</span>
      <span className="red-text">{props.errors.newPassword}</span>
      <h1>Reset Password</h1>
      <InputGroup className="mb-3">
        <InputGroupText>
          <CIcon icon={cilLockLocked} />
        </InputGroupText>
        <Input
          type="password"
          name="newpassword"
          placeholder="Password"
          onChange={(e) => {
            passwordReg.onChange(e) // method from hook form register
          }}
          innerRef={passwordReg.ref}
        />
      </InputGroup>
      {errors.newpassword && (
        <p className="edge-form-input-error">
          {errors.newpassword.message}
          <MyTooltip id="passwordRegister" tooltip={umTips.passwordHints} showTooltip={true} />
        </p>
      )}
      <InputGroup className="mb-4">
        <InputGroupText>
          <CIcon icon={cilLockLocked} />
        </InputGroupText>
        <Input
          type="password"
          name="newpassword2"
          placeholder="Repeat password"
          onChange={(e) => {
            password2Reg.onChange(e) // method from hook form register
          }}
          innerRef={password2Reg.ref}
        />
      </InputGroup>
      {errors.newpassword2 && (
        <p className="edge-form-input-error">{errors.newpassword2.message}</p>
      )}
      <Button color="primary" type="submit" className="px-4" block>
        Reset Password
      </Button>
    </Form>
  )
}

export default ResetPasswordForm
