#!/bin/bash
read -p "compress 0.0?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 0.0.zip 0.0*
fi
read -p "compress 0.1?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 0.1.zip 0.1*
fi
read -p "compress 0.2?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 0.2.zip 0.2*
fi
read -p "compress 0.3?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 0.3.zip 0.3*
fi
read -p "compress 0.4?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 0.4.zip 0.4*
fi
read -p "compress 0.5?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 0.5.zip 0.5*
fi
read -p "compress 0.6?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 0.6.zip 0.6*
fi
read -p "compress 0.7?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 0.7.zip 0.7*
fi
read -p "compress 0.8?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 0.8.zip 0.8*
fi
read -p "compress 0.9?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 0.9.zip 0.9*
fi
read -p "compress 1.0?" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip 1.0.zip 1.0*
fi

read -p "bundle into one archive" -n 1 -r
echo    # (optional) move to a new line
if [[ $REPLY =~ ^[Yy]$ ]]
then
    7z a -tzip -mx0 ../RIPPLE_1.zip *.zip
fi

#read -p "clean up?" -n 1 -r
#echo    # (optional) move to a new line
#if [[ $REPLY =~ ^[Yy]$ ]]
#then
    #rm *.zip
#fi
