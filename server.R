library(shiny)

shinyServer(function(input, output) {

  observe({
    result <- try(pow(rho1=ifelse(input$rho1_type==1, input$rho1, NA),
                      rho1_p=ifelse(input$rho1_type==1, NA, input$rho1),
                      rho2=ifelse(input$rho2_type==1, input$rho2, NA),
                      rho2_p=ifelse(input$rho2_type==1, NA, input$rho2),
                      delta=ifelse(input$delta_type==1, input$delta, NA),
                      delta_p=ifelse(input$delta_type==1, NA, input$delta),
                      alpha1=input$alpha1, alpha2=input$alpha2, n=input$n))
    
    if(class(result) == "try-error"){
      output$error <- renderUI({div(result, style = "color:red")})
      
      output$pBoth <- renderText("")
      output$pB1 <- renderText("")
      output$pB2 <- renderText("")
      output$pNone <- renderText("")
      
      output$rho1_simple <- renderText("")
      output$rho1_partial <- renderText("")
      output$rho2_simple <- renderText("")
      output$rho2_partial <- renderText("")
      output$delta_simple <- renderText("")
      output$delta_partial <- renderText("")
      
      output$scorr1 <- renderText("")
      output$pcorr1 <- renderText("")
      output$scorr2 <- renderText("")
      output$pcorr2 <- renderText("")
      output$r_x1x2 <- renderText("")
      output$pcorr3 <- renderText("")
    } else {
      output$error <- renderUI({})
      
      output$pBoth <- renderText(round(result["p2"], 3))
      output$pB1 <- renderText(round(result["p1.1"], 3))
      output$pB2 <- renderText(round(result["p1.2"], 3))
      output$pNone <- renderText(round(1 - result["p2"] - result["p1.1"] - result["p1.2"], 3))
      
      output$rho1_simple <- renderText(round(result["rho1"], 3))
      output$rho1_partial <- renderText(round(result["rho1_p"], 3))
      output$rho2_simple <- renderText(round(result["rho2"], 3))
      output$rho2_partial <- renderText(round(result["rho2_p"], 3))
      output$delta_simple <- renderText(round(result["delta"], 3))
      output$delta_partial <- renderText(round(result["delta_p"], 3))
      
      output$scorr1 <- renderText(round(result["scorr1"], 3))
      output$pcorr1 <- renderText(round(result["pcorr1"], 3))
      output$scorr2 <- renderText(round(result["scorr2"], 3))
      output$pcorr2 <- renderText(round(result["pcorr2"], 3))
      output$r_x1x2 <- renderText(round(result["r_x1x2"], 3))
      output$pcorr3 <- renderText(round(result["pcorr3"], 3))
    }
  }) # end observe
   
}) # end call to shinyServer()
