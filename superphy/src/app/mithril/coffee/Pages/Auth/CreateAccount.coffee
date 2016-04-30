# coffeelint: disable=max_line_length

class CreateAccount extends Page
    Routes.add('/SignUp', @)
    class FormGroup
        @controller: (args) ->
            args = args or {}
            @placeholder = args.placeholder or ""
            @type = args.type or "text"
            @value = args.value or m.prop("")
            @help_block = args.help_block or null
            return @
        @view: (ctrl) ->
            m('.form-group'
                m('label.col-sm-2 control-label', ctrl.placeholder)
                    m('.col-sm-8'
                        m('input.form-control', {
                            type:ctrl.type
                            oninput:m.withAttr('value', ctrl.value)
                            value:ctrl.value()
                            placeholder: ctrl.placeholder
                        })
                        if ctrl.help_block
                            m('span.help-block', "(#{ctrl.help_block})")
                    )
            )

    @controller: (args) ->
        @user = new User()
        @submit = () ->
            User.sign_up(@user)
            
        return @
    @view: (ctrl) ->
        super m('.container'
            m('.content-header', m('h1'
                m('span.title_part1', "CREATE"),
                m('span.title_part2', "ACCOUNT")
            ))
            m('.panel panel-default', m('panel-body', m('form.form-horizontal'
                m('.form-group', m('.col-sm-offset-2 col-sm-8'
                    m('p', "Already registered or want to make changes to your existing account?",
                        m('a', {tabindex:"9", href:"edit_account"}, "Sign in")
                    )
                    m('p', "To create a new account, fill in all fields below.")
                ))
                m.component(FormGroup, {
                    value: ctrl.user.username
                    type:"text"
                    placeholder:"Username"
                    help_block:"alphanumeric"
                })
                m.component(FormGroup, {
                    value: ctrl.user.password
                    type:"Password"
                    placeholder:"Password"
                    help_block:"6-10 characters"
                })
                m.component(FormGroup, {
                    value: ctrl.user.password2
                    type:"Password"
                    placeholder:"Re-enter Password"
                })
                m.component(FormGroup, {
                    value: ctrl.user.first_name
                    type:"text"
                    placeholder:"First Name"
                })
                m.component(FormGroup, {
                    value: ctrl.user.last_name
                    id:"inputLastname"
                    type:"text"
                    placeholder:"Last Name"
                })
                m.component(FormGroup, {
                    value: ctrl.user.email
                    id:"inputEmail"
                    type:"text"
                    placeholder:"Email"
                })
                m('.form-group'
                    m('.col-sm-offset-2 col-sm-8'
                        m("button.btn btn-primary[type=button]"
                            { onclick: -> ctrl.submit() }
                            "Create Account"
                        )
                    )
                )
            )))
        )