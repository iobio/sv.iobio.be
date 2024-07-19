# Use an official Node.js runtime as a parent image
FROM node:18

USER root

# Set the working directory
WORKDIR /

# Copy package.json and package-lock.json
COPY package*.json ./

# Install dependencies
# Install SQLite
RUN apt-get update
RUN apt-get install -y sqlite3 bcftools

#Node pagckage manager
RUN npm install

# Copy the rest of the application
COPY app.js ./app.js
COPY /data ./data
COPY /parsers ./parsers
COPY /testing ./testing

RUN chmod -R 777 /data
RUN chmod -R 777 /parsers
RUN chmod -R 777 /testing

# Expose the port the app runs on
EXPOSE 7477

# Command to run the app
CMD ["node", "app.js"]