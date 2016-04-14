# coffeelint: disable=max_line_length

class LoginForm extends Page
    Routes.add('/Login', @)
    @controller: (args) ->
        @user = new User()
        @submit = () ->
            User.log_in(@user)
        return @
    @view: (ctrl) ->
        super m('.container'
            m('.panel panel-default', m('panel-body', m('form.form-horizontal'
                m('h3', "Log in")
                m('input[type=text]', {
                    oninput:m.withAttr('value', ctrl.user.username)
                    value:ctrl.user.username()
                    placeholder:"Username"
                })
                m('input[type=Password]', {
                    oninput:m.withAttr('value', ctrl.user.password)
                    value:ctrl.user.password()
                    placeholder:"Password"
                })
                m("button.btn btn-primary[type=button]"
                    { onclick: -> ctrl.submit() }
                    "Login"
                )
            )))
        )