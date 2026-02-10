import React, { useState, useEffect } from 'react'
import { Col, Row } from 'reactstrap'
import { useForm } from 'react-hook-form'
import MySelect from '../../common/MySelect'
import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'

export const MultSelectInput = (props) => {
  const componentName = 'multSelectInput'
  const [form, setState] = useState({ ...components[componentName] })
  const [doValidation, setDoValidation] = useState(0)

  const {
    register,
    setValue,
    formState: { errors },
    trigger,
  } = useForm({
    mode: defaults['form_mode'],
  })

  useEffect(() => {
    form.selections = props.value ? props.value : []
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  const handleSelectChange = (selected) => {
    if (selected) {
      form.selections = selected
      if (form.selections.length > 0) {
        setValue('selections-hidden', 'valid', { shouldValidate: true })
      } else {
        setValue('selections-hidden', '', { shouldValidate: true })
      }
    } else {
      form.selections = []
      setValue('selections-hidden', '', { shouldValidate: true })
    }
    setDoValidation(doValidation + 1)
  }

  useEffect(() => {
    form.selections = props.value ? props.value : []
    if (form.selections.length > 0) {
      setValue('selections-hidden', 'valid', { shouldValidate: true })
    } else {
      setValue('selections-hidden', '', { shouldValidate: true })
    }
  }, []) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    //validate form
    trigger().then((result) => {
      form.validForm = result
      if (props.isOptional) {
        form.validForm = true
      }
      //force updating parent's inputParams
      props.setParams(form, props.name)
    })
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md="3">
          {props.tooltip ? (
            <MyTooltip
              id={`multSelctTooltip-${props.name}`}
              tooltip={props.tooltip}
              text={props.text}
              place={props.tooltipPlace ? props.tooltipPlace : defaults.tooltipPlace}
              color={props.tooltipColor ? props.tooltipColor : defaults.tooltipColor}
              showTooltip={props.showTooltip ? props.showTooltip : defaults.showTooltip}
              clickable={props.tooltipClickable ? props.tooltipClickable : false}
            />
          ) : (
            <>{props.text}</>
          )}
        </Col>
        <Col xs="12" md="9">
          <MySelect
            options={props.options}
            value={props.value}
            isMulti={true}
            placeholder={props.placeholder}
            checkbox={true}
            closeMenuOnSelect={false}
            hideSelectedOptions={false}
            onChange={handleSelectChange}
          />
          <input
            type="hidden"
            name={'selections-hidden'}
            id={'selections-hidden'}
            {...register('selections-hidden', {
              required: props.errMsg ? props.errMsg : 'At least one selection is required.',
            })}
          />
          {!props.isOptional && errors && errors['selections-hidden'] && (
            <p className="edge-form-input-error">{errors['selections-hidden'].message}</p>
          )}
        </Col>
      </Row>
      <br></br>
    </>
  )
}
