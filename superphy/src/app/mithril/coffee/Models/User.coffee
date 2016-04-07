# coffeelint: disable=max_line_length


User = (data) ->
    data = data or {}
    @username = m.prop(data.username or "")
    @password = m.prop(data.password or "")
    @password2 = m.prop(data.password2 or "")
    @first_name = m.prop(data.first_name or "")
    @last_name = m.prop(data.last_name or "")
    @email = m.prop(data.email or "")
    return
User.sign_up = (user) ->
    console.log("Sign Up Called #{JSON.stringify(user)}")
    #Do local error checking
    #Post to remote
    #Resolve response
User.log_in = (user) ->
    console.log("Log In Called #{JSON.stringify(user)}")
    #Do local error checking
    #Post to remote
    #Resolve Response
