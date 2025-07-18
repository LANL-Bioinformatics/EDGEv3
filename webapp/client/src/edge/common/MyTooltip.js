import React, { useState } from 'react'
import { Tooltip } from 'react-tooltip'
import { FaInfoCircle } from 'react-icons/fa'
import { HtmlText } from './HtmlText'

export const MyTooltip = (props) => {
  const [showTooltip, setShowTooltip] = useState(props.showTooltip)

  return (
    <>
      <a
        data-tooltip-id={props.id}
        className={props.className}
        onMouseOver={(e) => setShowTooltip(true)}
        onMouseLeave={(e) => {
          if (!props.showTooltip) setShowTooltip(false)
        }}
      >
        {props.text} &nbsp;
        <span style={{ color: props.color }} className={showTooltip ? '' : 'hide'}>
          <FaInfoCircle />
        </span>
      </a>
      <Tooltip
        style={{ zIndex: 9999 }}
        id={props.id}
        type={props.type}
        place={props.place ? props.place : 'right'}
        clickable={props.clickable ? props.clickable : false}
        offset={10}
      >
        <div
          style={{
            maxWidth: '600px',
          }}
        >
          <HtmlText text={props.tooltip} />
        </div>
      </Tooltip>
    </>
  )
}

export const WarningTooltip = (props) => {
  return (
    <>
      <a data-tooltip-id={'warning-' + props.id}>
        <span className="edge-form-input-error"></span>
      </a>
      <Tooltip
        wrapper="span"
        id={'warning-' + props.id}
        type={props.type ? props.type : 'warning'}
        place={props.place ? props.place : 'right'}
        className="edge-tooltip"
      >
        <HtmlText text={props.tooltip} />
      </Tooltip>
    </>
  )
}

export const ErrorTooltip = (props) => {
  return (
    <>
      <a data-tooltip-id={'error-' + props.id}>
        <span className="edge-text-size-large red-text">*</span>
      </a>
      <Tooltip
        className="edge-tooltip"
        id={'error-' + props.id}
        type={'error'}
        place={props.place ? props.place : 'right'}
      >
        <HtmlText text={props.tooltip} />
      </Tooltip>
    </>
  )
}

export default MyTooltip
