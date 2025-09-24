import React from 'react'
import { Badge } from 'reactstrap'
import UserTable from './UserTable'

const Users = (props) => {
  return (
    <div className="animated fadeIn">
      <div className="clearfix">
        <Badge color="danger" pill>
          Admin tool
        </Badge>
      </div>
      <br></br>
      <UserTable title={'Manage Users'} {...props} />
      <br></br>
    </div>
  )
}

export default Users
