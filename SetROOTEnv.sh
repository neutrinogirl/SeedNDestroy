#!/bin/bash

echo "Setting up ROOT_INCLUDE_PATH"
if [ -z "$ROOT_INCLUDE_PATH" ]; then
  export ROOT_INCLUDE_PATH=${PWD}/include
else
  export ROOT_INCLUDE_PATH=${ROOT_INCLUDE_PATH}:${PWD}/include
fi
