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
    #Do local error checking
    #Post to remote
    URL = 'api/users'
    response = m.request(
        headers: get_headers()
        method: "POST"
        url: "http://#{location.hostname}:5000/#{URL}"
        data: {
            'username': user.username()
            'password': user.password()
        }
        datatype: 'json'
        type: (response) ->
            return response
    )

User.log_in = (user) ->
    user = user or new User
    URL = 'api/token'
    auth = "Basic " + btoa("#{user.username()}:#{user.password()}")
    response = m.request(
        method: "GET"
        url: "http://#{location.hostname}:5000/#{URL}"
        config: (xhr) ->
            xhr.setRequestHeader('Content-Type', 'application/json')
            xhr.setRequestHeader('Authorization', auth)
 
        type: (response) ->
            if typeof(response) is "string"
                throw response
            else
                return response
    )