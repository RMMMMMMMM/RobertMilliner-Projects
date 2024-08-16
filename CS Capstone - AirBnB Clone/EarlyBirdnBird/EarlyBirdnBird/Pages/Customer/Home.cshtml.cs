using DataAccess;
using Infrastructure.Models;
using Microsoft.AspNetCore.Identity;
using Microsoft.AspNetCore.Mvc;
using Microsoft.AspNetCore.Mvc.RazorPages;
using Microsoft.EntityFrameworkCore;
using System.ComponentModel.DataAnnotations;

namespace EarlyBirdnBird.Pages.Customer
{
    public class HomeModel : PageModel
    {
        private readonly ApplicationDbContext _dbContext;
        private readonly UserManager<AppUser> _userManager;

        public HomeModel(ApplicationDbContext db, UserManager<AppUser> userManager)
        {
            _dbContext = db;
            _userManager = userManager;
            Reservations = new List<Reservation>();
        }

        public IList<Reservation> Reservations { get; set; }

        public async Task OnGetAsync()
        {
            if (!User.Identity.IsAuthenticated)
            {
                RedirectToPage("/Account/Login", new { area = "Identity" });
                return;
            }

            var currentUser = User.FindFirst(System.Security.Claims.ClaimTypes.NameIdentifier)?.Value;

            Reservations = _dbContext.Reservation
                .Where(x => x.AppUser.Id == currentUser)
                .Select(Reservation => new Reservation
                {
                    Id = Reservation.Id,
                    Rental = Reservation.Rental,
                    IndividualTotal = Reservation.IndividualTotal,
                    SalesTax = Reservation.SalesTax,
                    GrandTotal = Reservation.GrandTotal,
                    Guests = Reservation.Guests,
                    Status = Reservation.Status,
                    BookedStart = Reservation.BookedStart,
                    BookedEnd = Reservation.BookedEnd,
                })
                .ToList();
        }

        public async Task<IActionResult> OnPostAsync(int id)
        {

            //NEED ERROR HANDLING!!!!!!!!!!!!!!!!!!!

            var currentUser = await _userManager.GetUserAsync(User);
            var currentId = User.FindFirst(System.Security.Claims.ClaimTypes.NameIdentifier)?.Value;

            if (currentUser == null)
            {
                // User not found, handle the error accordingly
                return NotFound("User not found.");
            }

            //find reservation and rental/host

            var rentalId = _dbContext.Reservation
                .Where(x => x.Id == id)
                .Select(s => s.Rental)
                .SingleOrDefault();
            var rentalSD = _dbContext.Reservation
                .Where(x => x.Id == id)
                .Select(s => s.BookedStart)
                .SingleOrDefault();
            var rentalED = _dbContext.Reservation
                .Where(x => x.Id == id)
                .Select(s => s.BookedEnd)
                .SingleOrDefault();
            var rentalName = _dbContext.Rentals
                .Where(x => x.Id == rentalId)
                .Select(s => s.Name)
                .SingleOrDefault();

            var rentalOwner = _dbContext.Rentals
                .Where(x => x.Id == rentalId)
                .Select(s => s.OwnerName)
                .SingleOrDefault();

            var Rental = _dbContext.Rentals.Where(r => r.Id == rentalId).Single();

            
            var reservation = await _dbContext.Reservation.FindAsync(id);
            // Update the status of the reservation to "3" (Cancelled)
            reservation.Status = 3;

            // Save the changes to the database
            await _dbContext.SaveChangesAsync();

            string custsubject = "Reservation Cancellation";
            string custbody = "Reservation for " + rentalName + " has been Cancelled on these dates: " + rentalSD + " - " + rentalED + "";

            UserMail custMail = new UserMail()
            {
                fromUser = currentUser,
                toUser = rentalOwner.Email,
                subject = custsubject,
                body = custbody,
                DateSent = DateTime.Now
            };

            _dbContext.UserMails.Add(custMail);
            await _dbContext.SaveChangesAsync();


            string hostsubject = "Reservation Cancellation";
            string hostbody = "Reservation for " + rentalName + " has been Cancelled on these dates: " + rentalSD + " - " + rentalED + ""
                + " You have been refunded: $" + reservation.GrandTotal;

            UserMail hostMail = new UserMail()
            {
                fromUser = Rental.OwnerName,
                toUser = currentUser.Email,
                subject = hostsubject,
                body = hostbody,
                DateSent = DateTime.Now
            };

            _dbContext.UserMails.Add(hostMail);
            await _dbContext.SaveChangesAsync();



            // Return to home
            return RedirectToPage("./Home");
        }
    }
}
