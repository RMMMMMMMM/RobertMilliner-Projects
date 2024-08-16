using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.EntityFrameworkCore;


namespace EarlyBirdnBird.Pages
{
    public class InboxModel : PageModel
    {
        private readonly ApplicationDbContext _dbContext;
        private readonly UserManager<AppUser> _userManager;

        public InboxModel(ApplicationDbContext dbContext, UserManager<AppUser> userManager)
        {
            _dbContext = dbContext;
            _userManager = userManager;
            UserMail = new List<UserMail>();
        }

        public IList<UserMail> UserMail { get; set; }

        //public async Task OnGetAsync()
        //{
        //    var currentUser = _userManager.GetUserName(User);
        //    var mail = from m in _dbContext.UserMails select m;
        //    mail = mail.Where(x => x.toUser.Equals(currentUser, StringComparison.OrdinalIgnoreCase));
        //    UserMail = await mail.ToListAsync();
        //}

        public async Task OnGetAsync()
        {
            if (!User.Identity.IsAuthenticated)
            {
                RedirectToPage("/Account/Login", new { area = "Identity" });
                return;
            }

            var currentUser = await _userManager.GetUserAsync(User);

            UserMail = _dbContext.UserMails
                .Where(x => x.toUser == currentUser.UserName)
                .Select(mail => new UserMail
                {
                    Id = mail.Id,
                    fromUser = mail.fromUser,
                    toUser = mail.toUser,
                    subject = mail.subject,
                    body = mail.body,
                    DateSent = mail.DateSent,
                    isRead = mail.isRead,
                    isDeleted = mail.isDeleted,
                })
                .ToList();
        }
    }
}
