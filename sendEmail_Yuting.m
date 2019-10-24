
function sendEmail_Yuting(subject)

props = java.lang.System.getProperties; 
props.setProperty( 'mail.smtp.starttls.enable', 'true' );
props.setProperty( 'mail.smtp.socketFactory.port', '587' );
props.setProperty('mail.smtp.auth','true'); 
setpref('Internet','SMTP_Server','smtp-mail.outlook.com'); 
setpref('Internet','SMTP_Username','yutingchen0604@hotmail.com'); 
setpref('Internet','SMTP_Password','AaBb14207'); 
% transport.connect("smtp.live.com",587,null,null);

sendmail('yutingchen0604@hotmail.com',subject)


end