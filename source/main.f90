PROGRAM MAIN
  USE DRIVERS
  USE KEYS, ONLY : ACTION, REPLICATESTRUCT, NBENDREPLIC, REPLICFILE, ENERGYRECALC
  USE HELIXUTILS, ONLY : OUTPUTREPLIC
  USE CHAINUTILS, ONLY : CHAIN
  
  IMPLICIT NONE
  TYPE(CHAIN), TARGET :: WLC
  TYPE(CHAIN), POINTER :: CHAINP

  CHAINP=> WLC

  ! read in keyword parameters
  CALL READKEY   

  ! run the appropriate driver subroutines
  SELECT CASE(ACTION)
  CASE('GETSTRUCT')     
     CALL SETUPDRIVER(CHAINP)
     CALL GETSTRUCTDRIVER(CHAINP)
     IF (REPLICATESTRUCT) CALL OUTPUTREPLIC(CHAINP,NBENDREPLIC,REPLICFILE)
     CALL CLEANUPDRIVER(CHAINP)
  CASE('OPTIMIZE')
     CALL SETUPDRIVER(CHAINP)
     CALL OPTIMIZEDRIVER(CHAINP)
     IF (REPLICATESTRUCT) CALL OUTPUTREPLIC(CHAINP,NBENDREPLIC,REPLICFILE)
     CALL CLEANUPDRIVER(CHAINP)
  CASE('BASINHOP')
     CALL SETUPDRIVER(CHAINP)
     CALL BASINHOPDRIVER(CHAINP)
     CALL CLEANUPDRIVER(CHAINP)
  CASE('DATABASEPARSE')
     CALL DBPARSEDRIVER
  CASE DEFAULT
     PRINT*, 'No known action specified. Doing nothing'
     PRINT*, 'ACTION: ', ACTION  
  END SELECT

END PROGRAM MAIN
