# Getting Server Running on CHPC

---

## Build the Container

`docker build --platform linux/amd64 -t emersonlebleu/sv_backend_server:amd64v2.2 .`

-   `--platform` option allows you to specify the platform to build for because CHPC runs on an x86 (AMD24) and most of our new machines are ARM64
-   `-t` allows you to name and tag your build, the content after the `:` is the tag

## Push the container to docker hub THEN --> Pull the Container

`singularity pull --name sv_be_v2.2.sif docker://emersonlebleu/sv_backend_server:amd64v2.2`

-   This command will create a `.sif` file within your current directory
-   Note the name may be not exactly the same as the docker name, by convention there aren't colons in the sif name

## Execute the `node /app/app.js`

`singularity exec sv_backend_server_amd64_2.2.sif node /app.js`

-   Once there is a sif file we can exec commands
-   Just running the sif should work but for some reason I couldn't get that working the file directories were wrong and it couldn't find the app.js
