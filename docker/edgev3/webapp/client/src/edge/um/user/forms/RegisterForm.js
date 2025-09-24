import React from 'react'
import { Form, Button, InputGroup, InputGroupText, Input } from 'reactstrap'
import { useForm } from 'react-hook-form'
import CIcon from '@coreui/icons-react'
import { cilUser, cilLockLocked } from '@coreui/icons'
import { MyTooltip } from '../../../common/MyTooltip'
import { umTips } from '../../defaults'

const RegisterForm = (props) => {
  const {
    register,
    handleSubmit,
    formState: { errors },
    watch,
  } = useForm({
    mode: 'onChange',
  })

  const firstNameReg = {
    ...register('firstName', {
      required: 'Please enter your first name',
      pattern: {
        value: /^[A-Za-z]+$/,
        message: 'Your name must be composed only with letters',
      },
    }),
  }

  const lastNameReg = {
    ...register('lastName', {
      required: 'Please enter your last name',
      pattern: {
        value: /^[A-Za-z]+$/,
        message: 'Your name must be composed only with letters',
      },
    }),
  }

  const emailReg = {
    ...register('email', {
      required: 'Email is required',
      pattern: {
        // Validation pattern
        value: /^[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,4}$/i,
        message: 'Invalid email address',
      },
    }),
  }

  const passwordReg = {
    ...register('password', {
      required: 'Please enter a password',
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
    ...register('confirmPassword', {
      validate: (value) => value === watch('password') || 'The passwords do not match',
    }),
  }

  return (
    <Form onSubmit={handleSubmit(props.onSubmit)}>
      <h1>Sign Up</h1>
      <p className="text-muted">Create your account</p>
      <InputGroup className="mb-3">
        <InputGroupText className="text-muted"> First Name </InputGroupText>
        <Input
          type="text"
          name="firstName"
          onChange={(e) => {
            firstNameReg.onChange(e) // method from hook form register
          }}
          innerRef={firstNameReg.ref}
        />
      </InputGroup>
      {errors.firstName && <p className="edge-form-input-error">{errors.firstName.message}</p>}
      <InputGroup className="mb-3">
        <InputGroupText className="text-muted"> Last Name </InputGroupText>
        <Input
          type="text"
          name="lastName"
          onChange={(e) => {
            lastNameReg.onChange(e) // method from hook form register
          }}
          innerRef={lastNameReg.ref}
        />
      </InputGroup>
      {errors.lastName && <p className="edge-form-input-error">{errors.lastName.message}</p>}
      <span className="red-text">{props.errors.email}</span>
      <InputGroup className="mb-3">
        <InputGroupText>
          <CIcon icon={cilUser} />
        </InputGroupText>
        <Input
          type="text"
          name="email"
          placeholder="Email"
          onChange={(e) => {
            emailReg.onChange(e) // method from hook form register
          }}
          innerRef={emailReg.ref}
        />
      </InputGroup>
      {errors.email && <p className="edge-form-input-error">{errors.email.message}</p>}

      <InputGroup className="mb-3">
        <InputGroupText>
          <CIcon icon={cilLockLocked} />
        </InputGroupText>
        <Input
          type="password"
          name="password"
          placeholder="Password"
          onChange={(e) => {
            passwordReg.onChange(e) // method from hook form register
          }}
          innerRef={passwordReg.ref}
        />
      </InputGroup>
      {errors.password && (
        <p className="edge-form-input-error">
          {errors.password.message}
          <MyTooltip id="passwordRegister" tooltip={umTips.passwordHints} showTooltip={true} />
        </p>
      )}
      <InputGroup className="mb-4">
        <InputGroupText>
          <CIcon icon={cilLockLocked} />
        </InputGroupText>
        <Input
          type="password"
          name="confirmPassword"
          placeholder="Repeat password"
          onChange={(e) => {
            password2Reg.onChange(e) // method from hook form register
          }}
          innerRef={password2Reg.ref}
        />
      </InputGroup>
      {errors.confirmPassword && (
        <p className="edge-form-input-error">{errors.confirmPassword.message}</p>
      )}
      <Button color="primary" disabled={props.loading} type="submit" block>
        Create Account
      </Button>
    </Form>
  )
}

export default RegisterForm
