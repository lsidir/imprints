list(APPEND experiments
    main
    )

foreach(ee ${experiments})
    add_executable(${ee} "${ee}.c")
    target_link_libraries(${ee} impsjoin-core )
endforeach()

