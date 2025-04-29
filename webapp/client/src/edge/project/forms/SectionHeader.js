import React, { useState } from 'react'
import { Button, ButtonGroup, CardHeader } from 'reactstrap'
import { cilChevronBottom, cilChevronTop } from '@coreui/icons'
import CIcon from '@coreui/icons-react'
import { colors, defaults } from '../../common/util'
import { ErrorTooltip } from '../../common/MyTooltip'

export const Header = (props) => {
  const [headerColor, setHeaderColor] = useState(colors.light)
  return (
    <CardHeader
      key={props.id}
      style={{ backgroundColor: headerColor, border: '1px' }}
      onMouseEnter={(e) => setHeaderColor(colors.secondary)}
      onMouseLeave={() => setHeaderColor(colors.light)}
    >
      {props.toggle ? (
        <>
          <span
            className="edge-card-header-action"
            data-target="#collapseParameters"
            onClick={props.toggleParms}
          >
            {props.collapseParms ? (
              <CIcon icon={cilChevronBottom} className="me-2" />
            ) : (
              <CIcon icon={cilChevronTop} className="me-2" />
            )}
          </span>
          &nbsp;&nbsp;&nbsp;
        </>
      ) : (
        <>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;</>
      )}
      <span className="edge-card-header-title">{props.title}</span>
      {!props.isValid && (
        <>
          &nbsp;
          <ErrorTooltip id={props.id} tooltip={props.errMessage} place={defaults.tooltipPlace} />
        </>
      )}
      {props.onoff && (
        <ButtonGroup style={{ float: 'right' }} className="mr-3" aria-label="First group" size="sm">
          <Button
            color="outline-primary"
            onClick={() => {
              props.setOnoff(true)
            }}
            active={props.paramsOn}
          >
            On
          </Button>
          <Button
            color="outline-primary"
            onClick={() => {
              props.setOnoff(false)
            }}
            active={!props.paramsOn}
          >
            Off
          </Button>
        </ButtonGroup>
      )}
    </CardHeader>
  )
}
