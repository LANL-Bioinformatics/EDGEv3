import React, { useState, useEffect } from 'react'
import { Button, ButtonGroup, Col, Row } from 'reactstrap'
import { MyTooltip } from '../../common/MyTooltip'
import { defaults } from '../../common/util'
import { components } from './defaults'

export const Switcher = (props) => {
  const componentName = 'switcher'
  const [form, setState] = useState({ ...components[componentName], isTrue: props.defaultValue })
  const [doValidation, setDoValidation] = useState(0)

  const setNewState2 = (name, value) => {
    setState({
      ...form,
      [name]: value,
    })
    setDoValidation(doValidation + 1)
  }

  useEffect(() => {
    setState({ ...components[componentName], isTrue: props.defaultValue })
  }, [props.reset]) // eslint-disable-line react-hooks/exhaustive-deps

  useEffect(() => {
    //force updating parent's inputParams
    props.setParams(form, props.name)
  }, [doValidation]) // eslint-disable-line react-hooks/exhaustive-deps

  return (
    <>
      <Row>
        <Col md="3">
          {props.tooltip ? (
            <MyTooltip
              id={`switcherTooltip-${props.name}`}
              tooltip={props.tooltip}
              text={props.text}
              place={props.tooltipPlace ? props.tooltipPlace : defaults.tooltipPlace}
              color={props.tooltipColor ? props.tooltipColor : defaults.tooltipColor}
              showTooltip={props.showTooltip ? props.showTooltip : defaults.showTooltip}
            />
          ) : (
            <>{props.text}</>
          )}{' '}
        </Col>
        <Col xs="12" md="9">
          <ButtonGroup className="mr-3" aria-label="First group" size="sm">
            <Button
              color="outline-primary"
              onClick={() => {
                setNewState2('isTrue', true)
              }}
              active={form.isTrue}
            >
              {props.trueText}
            </Button>
            <Button
              color="outline-primary"
              onClick={() => {
                setNewState2('isTrue', false)
              }}
              active={!form.isTrue}
            >
              {props.falseText}
            </Button>
          </ButtonGroup>
        </Col>
      </Row>
    </>
  )
}
