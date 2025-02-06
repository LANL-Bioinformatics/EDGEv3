import React, { useState, useEffect } from 'react'
import Select from 'react-select'
import { useSelector } from 'react-redux'
import { Button, Modal, ModalHeader, ModalBody, ModalFooter } from 'reactstrap'

import { getData, apis } from '../../common/util'
import startCase from 'lodash.startcase'

const UserSelector = (props) => {
  const user = useSelector((state) => state.user)
  const [options, setOptions] = useState([])

  //loaded only once
  useEffect(() => {
    const getUsers = () => {
      getData(apis.userUsers)
        .then((data) => {
          return data.users.reduce((options, obj) => {
            if (props.type === 'admin') {
              options.push({
                value: obj.email,
                label: obj.firstName + ' ' + obj.lastName + ', ' + obj.email,
                id: obj.id,
              })
            } else {
              if (obj.email !== user.profile.email) {
                options.push({
                  value: obj.email,
                  label: obj.firstName + ' ' + obj.lastName + ', ' + obj.email,
                  id: obj.id,
                })
              }
            }
            return options
          }, [])
        })
        .then((options) => {
          setOptions(options)
        })
        .catch((err) => {
          alert(err)
        })
    }
    getUsers()
  }, [props, user])

  return (
    <Modal isOpen={props.isOpen} centered>
      <ModalHeader>Choose users to {props.action} the selected projects</ModalHeader>
      <ModalBody>
        <Select
          options={options}
          onChange={props.onChange}
          isMulti={true}
          centered={true}
          autoFocus={true}
        />
      </ModalBody>
      <ModalFooter>
        <Button color="primary" onClick={props.handleClickYes}>
          {startCase(props.action)}
        </Button>{' '}
        <Button color="secondary" onClick={props.handleClickClose}>
          Cancel
        </Button>
      </ModalFooter>
    </Modal>
  )
}

export default UserSelector
