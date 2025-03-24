import React, { useState, useEffect } from 'react'
import { Col, Row, Input } from 'reactstrap'
import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'

export const FileUpload = (props) => {
  const componentName = 'fileUpload'
  const inputStyle = defaults.inputStyle
  const inputStyleWarning = defaults.inputStyleWarning
  const [form, setState] = useState({ ...components[componentName] })
  const [doValidation, setDoValidation] = useState(0)

  useEffect(() => {
    //validate form
    if (form.file || props.isOptional) {
      form.errMessage = ''
      form.validForm = true
    } else {
      form.errMessage = 'File is required.'
      form.validForm = false
    }
    //force updating parent's inputParams
    props.setParams(form, props.name)
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md="3">
          {props.tooltip ? (
            <MyTooltip
              id={`fileUploadTooltip-${props.name}`}
              tooltip={props.tooltip}
              text={props.text}
              place={props.tooltipPlace ? props.tooltipPlace : defaults.tooltipPlace}
              color={props.tooltipColor ? props.tooltipColor : defaults.tooltipColor}
              showTooltip={props.showTooltip ? props.showTooltip : defaults.showTooltip}
            />
          ) : (
            <>{props.text}</>
          )}
        </Col>
        <Col xs="12" md="9">
          <Input
            style={
              // eslint-disable-next-line prettier/prettier
              (form.validForm || props.isValidFileInput)
                ? inputStyle
                : inputStyleWarning
            }
            type="file"
            id="file"
            onChange={(event) => {
              form['file'] = event.target.files[0]
              setDoValidation(doValidation + 1)
            }}
            accept={props.accept ? props.accept : '*'}
          />
        </Col>
      </Row>
    </>
  )
}
