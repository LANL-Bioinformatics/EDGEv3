import React from 'react'
import { Form, Input, Button, InputGroup, InputGroupText } from 'reactstrap'
import { useForm } from 'react-hook-form'
import CIcon from '@coreui/icons-react'
import { cilUser } from '@coreui/icons'

const ResetPasswordLinkForm = (props) => {
  const {
    register,
    handleSubmit,
    formState: { errors },
  } = useForm({
    mode: 'onChange',
  })

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

  return (
    <Form onSubmit={handleSubmit(props.onSubmit)}>
      <span className="red-text">{props.messages.getResetPasswordLink}</span>
      <span className="red-text">{props.errors.getResetPasswordLink}</span>
      <h1>Forgot Password</h1>
      <div>
        <span className="red-text">{props.errors.resetPasswordLink}</span>
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
        <Button color="primary" type="submit" className="px-4" block>
          Get ResetPassword Link
        </Button>
      </div>
    </Form>
  )
}

export default ResetPasswordLinkForm
