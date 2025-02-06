## INSTALLATION PREREQUISITES

### Install Node20
https://nodejs.org/dist/latest-v20.x/

### Install pm2
`npm install pm2@latest -g`

### Install MongoDB Community Edition
https://docs.mongodb.com/manual/installation/#mongodb-community-edition-installation-tutorials

## INSTALL webapp

1. Move/copy EDGEv3 folder to the installation directory

2. Inside EDGEv3/installation folder, run the installation script 

    `./install.sh`

3. Create environment variables

    The web client and web server each rely on environment variables for their configuration.
    You can define those environment variables in `.env` files.

    Here's how you can define them in `.env` files:

    - Populate the "client build" environment configuration file (i.e. `webapp/client/.env`). 

        You can initialize it based upon the corresponding example file:
        ```shell
        cp webapp/client/.env.example \
        webapp/client/.env
        ```
        > Those environment variables are used within `webapp/client/src/config.js`.

    -  Create a build directory with a production build of the client:
        ```shell
        cd webapp/client
        npm run build
        ```

    - Populate the server environment configuration file (i.e. `webapp/server/.env`). 
    
        You can initialize it based upon the corresponding example file:
        ```shell
        cp webapp/server/.env.example \
        webapp/server/.env
        ```
        > Those environment variables are used within `webapp/server/config.js`.

## START webapp

1. Start MongoDB if it's not started yet

2. Inside EDGEv3 folder, run the pm2 start command 

    `pm2 start pm2.config.js`
    
## STOP webapp

    pm2 stop all
