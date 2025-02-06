import React from 'react'
import { CFooter } from '@coreui/react'

import lanlLogo from 'src/assets/images/logo_LANL.png'
import doeLogo from 'src/assets/images/DOE-logo.png'
import nnsaLogo from 'src/assets/images/logo_NNSA.png'
import dtraLogo from 'src/assets/images/logo_DTRA.png'

const AppFooter = () => {
  return (
    <CFooter className="edge-footer">
      <div>
        <span className="edge-footer-item">
          <a target="_lanl" href="http://www.lanl.gov">
            <img alt="lanl logo" src={lanlLogo} style={{ width: 150, height: 30 }} />
          </a>
        </span>
      </div>
      <span className="edge-footer-item">
        <div className="edge-text-center edge-text-size-small">
          <a target="_lanl" href="https://www.lanl.gov//resources/web-policies/copyright-legal.php">
            <small>Terms of Use, Privacy</small>
          </a>
        </div>
        <p className="edge-text-center edge-text-size-small">
          <small>
            Managed by Triad National Security, LLC for the U.S Dept. of Energy&apos;s NNSA
          </small>
          <br></br>
          <small>© Copyright Triad National Security, LLC. All Rights Reserved.</small>
        </p>
      </span>
      <div className="edge-footer-item">
        <span className="mr-1">
          <a target="_new" href="http://www.dtra.mil">
            <img alt="dtra logo" style={{ width: 50, height: 50 }} src={dtraLogo} />
          </a>
          <a target="_new" href="https://www.energy.gov/">
            <img alt="doe logo" style={{ width: 50, height: 50 }} src={doeLogo} />
          </a>
          <a target="_new" href="http://nnsa.energy.gov/">
            <img alt="nnsa logo" style={{ width: 70, height: 50 }} src={nnsaLogo} />
          </a>
          {/* <a target="_new" href="https://access-ci.org/">
            <img
              alt="access logo"
              style={{ width: 170, height: 50 }}
              src={'https://access-ci.org/wp-content/uploads/2022/08/access-logo-footer.svg'}
            />
          </a> */}
        </span>
      </div>
    </CFooter>
  )
}

export default React.memo(AppFooter)
